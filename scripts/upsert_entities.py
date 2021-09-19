import argparse
import pandas as pd
import requests

from oauth2client.client import GoogleCredentials


# function to get authorization bearer token for requests
def get_access_token():
    """Get access token."""

    scopes = ["https://www.googleapis.com/auth/userinfo.profile", "https://www.googleapis.com/auth/userinfo.email"]
    credentials = GoogleCredentials.get_application_default()
    credentials = credentials.create_scoped(scopes)

    return credentials.get_access_token().access_token


def call_rawls_batch_upsert(workspace_name, project, request):
    """Post entities to Terra workspace using batchUpsert."""

    # rawls request URL for batchUpsert
    uri = f"https://rawls.dsde-prod.broadinstitute.org/api/workspaces/{project}/{workspace_name}/entities/batchUpsert"

    # Get access token and and add to headers for requests.
    # -H  "accept: */*" -H  "Authorization: Bearer [token] -H "Content-Type: application/json"
    headers = {"Authorization": "Bearer " + get_access_token(), "accept": "*/*", "Content-Type": "application/json"}

    # capture response from API and parse out status code
    response = requests.post(uri, headers=headers, data=request)
    status_code = response.status_code

    if status_code != 204:  # entities upsert fail
        print(f"WARNING: Failed to upsert entities to {workspace_name}.")
        print(response.text)
        return

    # entities upsert success
    print(f"Successful upsert of entities to {workspace_name}.")


def write_request_json(request):
    """Create output file with json request."""

    save_name = "batch_upsert_request.json"
    with open(save_name, "w") as f:
        f.write(request)


def create_upsert_request(tsv):
    """Generate the request body for batchUpsert API."""

    df_tsv = pd.read_csv(tsv, sep="\t")

    # check if the tsv is formatted correctly - exit script if not in right load format
    entity_type_col_name = df_tsv.columns[0]
    if not entity_type_col_name.startswith("entity:"):
        print("Not a valid tsv. The .tsv does not start with column entity:[table_name]_id. Please correct and try again.")
        return

    # define which columns are array values and which are not
    single_attr_cols = ["sample_id", "recalibrated_bam", "recalibrated_bam_index",
                       "table_before", "table_after", "plots"]

    # templates for request body components
    template_req_body = '''[{"name":"VAR_ENTITY_ID",
                             "entityType":"VAR_ENTITY_TYPE",
                             "operations":[OPERATIONS_LIST]}]'''

    template_make_list_attr = '{"op":"CreateAttributeValueList","attributeName":"VAR_ATTRIBUTE_LIST_NAME"},'
    template_add_list_member = '{"op":"AddListMember","attributeListName":"VAR_ATTRIBUTE_LIST_NAME", "newMember":"VAR_LIST_MEMBER"},'
    template_make_single_attr = '{"op":"AddUpdateAttribute","attributeName":"VAR_ATTRIBUTE_NAME", "addUpdateAttribute":"VAR_ATTRIBUTE_MEMBER"},'

    # initiate string to capture all operation requests
    all_operation_requests = ''''''

    # for each column (none are arrays)
    for col in single_attr_cols:
        # get valuue in col from df
        attr_val = str(df_tsv.iloc[0][col])

        # add the request for attribute to list
        all_operation_requests += template_make_single_attr.replace("VAR_ATTRIBUTE_MEMBER", attr_val).replace("VAR_ATTRIBUTE_NAME", col)

    # remove trailing comma from the last request template
    all_operation_requests = all_operation_requests[:-1]

    # get the entity_type (table name) and entity_id (row id - must be unique)
    entity_id = str(df_tsv.iloc[0][0])
    entity_type = entity_type_col_name.rsplit("_", 1)[0].split(":")[1]

    # put entity_type and entity_id in request body template
    final_request = template_req_body.replace("VAR_ENTITY_ID", entity_id).replace("VAR_ENTITY_TYPE", entity_type)
    # put operations list into request body template
    final_request = final_request.replace("OPERATIONS_LIST", all_operation_requests)

    # write out a json of the request body
    write_request_json(final_request)

    return final_request


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-w', '--workspace_name', required=True, help='name of workspace in which to make changes')
    parser.add_argument('-p', '--project', required=True, help='billing project (namespace) of workspace in which to make changes')
    parser.add_argument('-t', '--tsv', required=True, help='.tsv file formatted in load format to Terra UI')
    args = parser.parse_args()

    # create request body for batchUpsert
    request = create_upsert_request(args.tsv)
    # call batchUpsert API (rawls)
    call_rawls_batch_upsert(args.workspace_name, args.project, request)
